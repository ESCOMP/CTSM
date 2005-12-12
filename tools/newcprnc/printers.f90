module printers

   use prec,          only: r8
   use types,         only: varinfo, diffinfo
   use specialvalues, only: uninit_int
   use nfwrappers,    only: wrap_inq_varid, wrap_get_var_int, wrap_get_vara_text, wrap_get_vara_int
   use netcdfids,     only: dateid, datesecid, nstephid, date_writtenid, time_writtenid

   implicit none

   include 'netcdf.inc'

PRIVATE

   public :: prheader                 ! print header info
   public :: prhddiff                 ! print header diffs
   public :: print_stats              ! print variable stats
   public :: print_comparison_stats   ! print variable difference stats

CONTAINS

   subroutine prheader (nfid, numcases, t, ntime)
!-------------------------------------------------------------------------------------------
! Purpose: Print header info for 1 or 2 cases
!
! Method:  Make netcdf calls to retrieve the values or print the value of the input argument
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: nfid(2)       ! file id(s)
      integer, intent(in) :: numcases      ! number of files
      integer, intent(in) :: t(2)          ! time index
      integer, intent(in) :: ntime(2)      ! number of time samples
!
! Local workspace
!
      character(len=8) :: date_written(2)  ! date stamp on the time sample
      character(len=8) :: time_written(2)  ! time of day stamp

      integer :: date(2)                   ! simulation date
      integer :: datesec(2)                ! simulation time of day
      integer :: nsteph(2)                 ! time step
      integer :: startc(2)                 ! start array for nf_get_vara_text calls
      integer :: kountc(2)                 ! kount array for nf_get_vara_text calls
      integer :: start(1)                  ! start array for nf_get_vara_int calls
      integer :: kount(1)                  ! kount array for nf_get_vara_int calls
      integer :: n                         ! file index
!
! Initialize variables in case they do not exist on the file
!
      date(:)    = uninit_int
      datesec(:) = uninit_int
      nsteph(:)  = uninit_int
      date_written(:) = 'N/A'
      time_written(:) = 'N/A'

      write(6,*)'SUMMARY OF HEADER INFO TIME SAMPLES', t(1), t(2),':'
      do n=1,numcases
         startc(:) = (/1,t(n)/)
         kountc(:) = (/8,1/)
         if (date_writtenid(n) > 0) call wrap_get_vara_text (nfid(n), date_writtenid(n), startc, kountc, date_written(n))
         if (time_writtenid(n) > 0) call wrap_get_vara_text (nfid(n), time_writtenid(n), startc, kountc, time_written(n))

         start(1) = t(n)
         kount(1) = 1
         if (nstephid(n) > 0)  call wrap_get_vara_int (nfid(n), nstephid(n), start, kount, nsteph(n))
         if (dateid(n) > 0)    call wrap_get_vara_int (nfid(n), dateid(n), start, kount, date(n))
         if (datesecid(n) > 0) call wrap_get_vara_int (nfid(n), datesecid(n), start, kount, datesec(n))
      end do

      if (numcases > 1) then
         write(6,800)date_written(1), date_written(2), &
                     time_written(1), time_written(2), &
                     t(1), t(2), &
                     ntime(1), ntime(2), &
                     nsteph(1), nsteph(2), &
                     date(1), date(2), &
                     datesec(1),  datesec(2)
      else
         write(6,801)date_written(1), &
                     time_written(1), &
                     t(1), &
                     ntime(1), &
                     nsteph(1), &
                     date(1), &
                     datesec(1)
      end if
      write(6,*)

      return

800   format(' date_written:                ',2(A8,1X),/, &
             ' time_written:                ',2(A8,1X),/, &
             ' time sample:                 ',2I9,/,      &
             ' total possible time samples: ',2I9,/,      &
             ' NSTEPH:                      ',2I9,/,      &
             ' DATE:                        ',2I9,/,      &
             ' DATESEC:                     ',2I9)

801   format(' date_written:                ',1(A8,1X),/, &
             ' time_written:                ',1(A8,1X),/, &
             ' time sample:                 ',1I9,/,      &
             ' total possible time samples: ',1I9,/,      &
             ' NSTEPH:                      ',1I9,/,      &
             ' DATE:                        ',1I9,/,      &
             ' DATESEC:                     ',1I9)
   end subroutine prheader

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine prhddiff (nfid, numcases)
!-------------------------------------------------------------------------------------------
! Purpose: Print header difference info between 2 files
!
! Method:  Make netcdf calls to retrieve the values. Print difference info if they do not
!          match
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: nfid(2)        ! file id(s)
      integer, intent(in) :: numcases       ! number of files
!
! Local workspace
!
      integer :: n                          ! file index
      integer :: ntrmid(2)                  ! netcdf id for ntrm
      integer :: ntrnid(2)                  ! netcdf id for ntrn
      integer :: ntrkid(2)                  ! netcdf id for ntrk
      integer :: ndbaseid(2)                ! netcdf id for ndbase
      integer :: nsbaseid(2)                ! netcdf id for nsbase
      integer :: nbdateid(2)                ! netcdf id for nbdate
      integer :: nbsecid(2)                 ! netcdf id for nbsec
      integer :: mdtid(2)                   ! netcdf id for mdt

      integer :: ntrm(2)       = uninit_int ! spectral truncation parameter M
      integer :: ntrn(2)       = uninit_int ! spectral truncation parameter N
      integer :: ntrk(2)       = uninit_int ! spectral truncation parameter K
      integer :: ndbase(2)     = uninit_int ! base day
      integer :: nsbase(2)     = uninit_int ! base seconds of base day
      integer :: nbdate(2)     = uninit_int ! base date
      integer :: nbsec(2)      = uninit_int ! base seconds of base date
      integer :: mdt(2)        = uninit_int ! time step
!
! Loop through files, get variable info (if available), and print if diffs are found
!
      do n=1,numcases
         call wrap_inq_varid (nfid(n), 'ntrm', ntrmid(n), -1)
         call wrap_inq_varid (nfid(n), 'ntrn', ntrnid(n), -1)
         call wrap_inq_varid (nfid(n), 'ntrk', ntrkid(n), -1)
         call wrap_inq_varid (nfid(n), 'ndbase', ndbaseid(n), -1)
         call wrap_inq_varid (nfid(n), 'nsbase', nsbaseid(n), -1)
         call wrap_inq_varid (nfid(n), 'nbdate', nbdateid(n), -1)
         call wrap_inq_varid (nfid(n), 'nbsec', nbsecid(n), -1)
         call wrap_inq_varid (nfid(n), 'mdt', mdtid(n), -1)
!
! Cannot use optional arg to get_var_int or it complains about scalar/array problems
!
         if (ntrmid(n) > 0)   call wrap_get_var_int (nfid(n), ntrmid(n), ntrm(n))
         if (ntrnid(n) > 0)   call wrap_get_var_int (nfid(n), ntrnid(n), ntrn(n))
         if (ntrkid(n) > 0)   call wrap_get_var_int (nfid(n), ntrkid(n), ntrk(n))
         if (ndbaseid(n) > 0) call wrap_get_var_int (nfid(n), ndbaseid(n), ndbase(n))
         if (nsbaseid(n) > 0) call wrap_get_var_int (nfid(n), nsbaseid(n), nsbase(n))
         if (nbdateid(n) > 0) call wrap_get_var_int (nfid(n), nbdateid(n), nbdate(n))
         if (nbsecid(n) > 0)  call wrap_get_var_int (nfid(n), nbsecid(n), nbsec(n))
         if (mdtid(n) > 0)    call wrap_get_var_int (nfid(n), mdtid(n), mdt(n))
      end do

      write(6,*)'SUMMARY OF TIME-INDEPENDENT HEADER DIFFERENCES:'
      if (ntrm(1)   /= ntrm(2)  ) write(6,100)'ntrm  :', ntrm(1)  , ntrm(2)  
      if (ntrn(1)   /= ntrn(2)  ) write(6,100)'ntrn  :', ntrn(1)  , ntrn(2)  
      if (ntrk(1)   /= ntrk(2)  ) write(6,100)'ntrk  :', ntrk(1)  , ntrk(2)  
      if (ndbase(1) /= ndbase(2)) write(6,100)'ndbase:', ndbase(1), ndbase(2)
      if (nsbase(1) /= nsbase(2)) write(6,100)'nsbase:', nsbase(1), nsbase(2)
      if (nbdate(1) /= nbdate(2)) write(6,100)'nbdate:', nbdate(1), nbdate(2)
      if (nbsec(1)  /= nbsec(2) ) write(6,100)'nbsec :', nbsec(1) , nbsec(2) 
      if (mdt(1)    /= mdt(2)   ) write(6,100)'mdt   :', mdt(1)   , mdt(2)   

      write(6,300)
      return

100   format(a,2i9)
300   format(132('*'))
   end subroutine prhddiff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine print_stats (varstats)
!-------------------------------------------------------------------------------------------
! Purpose: Print various statistics (max, min, etc.) contained in input
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      type(varinfo), intent(in) :: varstats   ! derived type containing all the info

      write(6,809) varstats%dimnames(1), varstats%dimnames(2), varstats%dimnames(3)
      write(6,810) varstats%imax, varstats%jmax, varstats%kmax, &
                   varstats%imin, varstats%jmin, varstats%kmin
      write(6,803) varstats%name, varstats%npossible, varstats%nfill, varstats%arrmax, varstats%arrmin
      write(6,805) varstats%avg
      
      return

803   format(1x,a8,1x,i6,1x,i6,1p,2e23.15)
805   format(10x,'avg abs field values:  ',1pe23.15,/)
809   format(24x,'(',a5,',',a5,',',a5,')')
810   format(24x,'(',i5,',',i5,',',i5,')', t47,'(',i5,',',i5,',',i5,')')
   end subroutine print_stats

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine print_comparison_stats (varstats, diffstats)
!-------------------------------------------------------------------------------------------
! Purpose: Print various statistics (number of diffs, biggest diff, etc.) contained in input
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      type(varinfo), intent(in) :: varstats(2)  ! derived type for single-variable stats
      type(diffinfo), intent(in) :: diffstats   ! derived type for difference stats
!
! Local workspace
!
      real(r8) :: digbar      ! mean number of digits of agreement
      real(r8) :: digworst    ! number of digits of agreement for worst poin in the domain

      write(6,809)varstats(1)%dimnames(1), varstats(1)%dimnames(2), varstats(1)%dimnames(3)

      if (diffstats%ndiff > 0) then
         write(6,810)varstats(1)%imax, varstats(1)%jmax, varstats(1)%kmax, &
                     varstats(1)%imin, varstats(1)%jmin, varstats(1)%kmin, &
                     diffstats%imax,  diffstats%jmax,  diffstats%kmax, &
                     diffstats%irmax, diffstats%jrmax, diffstats%krmax
         write(6,803)diffstats%name, diffstats%ndiff, varstats(1)%nfill, varstats(1)%arrmax, varstats(1)%arrmin, &
                        diffstats%dmax, diffstats%dvals(1), diffstats%rdmax, diffstats%rdvals(1), &
                     varstats(1)%npossible, varstats(2)%nfill, varstats(2)%arrmax, varstats(2)%arrmin, &
                        diffstats%dvals(2), diffstats%rdvals(2)
         write(6,812)varstats(2)%imax, varstats(2)%jmax, varstats(2)%kmax, &
                     varstats(2)%imin, varstats(2)%jmin, varstats(2)%kmin
!
! Compute # digits accuracy for worst case & avg differences
!
         digbar = diffstats%rdlnbar / diffstats%ndiff
         digworst = log10 (1. / diffstats%rdmax)

         write(6,805) varstats(1)%avg, diffstats%rms, diffstats%rdbar, &
                      varstats(2)%avg, digbar, digworst
      else
         write(6,810)varstats(1)%imax, varstats(1)%jmax, varstats(1)%kmax, &
                     varstats(1)%imin, varstats(1)%jmin, varstats(1)%kmin
         write(6,814)diffstats%name, diffstats%ndiff, varstats(1)%nfill, varstats(1)%arrmax, varstats(1)%arrmin, &
                        varstats(1)%npossible, varstats(2)%nfill, varstats(2)%arrmax, varstats(2)%arrmin
         write(6,812)varstats(2)%imax, varstats(2)%jmax, varstats(2)%kmax, &
                     varstats(2)%imin, varstats(2)%jmin, varstats(2)%kmin
         write(6,815)varstats(1)%avg, varstats(2)%avg
      end if
!
! The following is strictly for test-model: Actual mass-weighting not yet implemented
!
      write(6,'(a,a8,1pe11.4,/)') ' RMS ', diffstats%name, diffstats%rms
      return
!
! formats
!
803   format(1x,a8,1x,i6,1x,i6,1p,2e23.15,e8.1,e23.15,e8.1,e23.15,/, &
             10x,     i6,1x,i6,   2e23.15,8x  ,e23.15,8x,  e23.15)
809   format(24x,'(',a5,',',a5,',',a5,')')
810   format(24x,'(',i5,',',i5,',',i5,')', t47,'(',i5,',',i5,',',i5,')', &
             t79,'(',i5,',',i5,',',i5,')',t109,'(',i5,',',i5,',',i5,')')
812   format(24x,'(',i5,',',i5,',',i5,')', t47,'(',i5,',',i5,',',i5,')')
814   format(1x,a8,1x,i6,1x,i6,1p2e23.15,/, &
             10x,     i6,1x,i6,2e23.15)
805   format(10x,'avg abs field values:  ',1pe23.15,3x,'rms diff:',e8.1, &
              3x,'avg rel diff(npos): ',e8.1,/, &
             10x,'                       ',  e23.15,23x, &
             'avg decimal digits(ndif): ',0p,f4.1,' worst: ',f4.1)
815   format(10x,'avg abs field values:  ',1p,e23.15,/, &
             10x,'                       ',   e23.15)
   end subroutine print_comparison_stats
end module printers
