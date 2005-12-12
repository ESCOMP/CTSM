module gatherstats

   use prec,          only: r4, r8
   use types,         only: varinfo, diffinfo
   use specialvalues, only: uninit_int, uninit_r8
   use cntlvars,      only: iprs, ipre, jprs, jpre, kprs, kpre
   use utils,         only: endrun

   implicit none

PRIVATE

   include 'netcdf.inc'

   public :: gather_stats              ! gather non-comparison stats (max, min, avg, nfill)
   public :: gather_comparison_stats   ! gather comparison stats (max diff, max relative diff, ...)

CONTAINS

   subroutine gather_stats (name, ndims, dimsiz1, dimsiz2, dimsiz3, &
                            arr, fillvalue, doprint, dimnames, varstats)
!-------------------------------------------------------------------------------------------
! Purpose: Gather max, min, avg, and associated position info for this variable
!
! Method:  Check only those points where _FillValue does not apply. Copy output into
!          varstats
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      character(len=*), intent(in) :: name                         ! variable name
      integer, intent(in) :: ndims                                 ! number of dimensions
      integer, intent(in) :: dimsiz1                               ! 1st non-time dimension size
      integer, intent(in) :: dimsiz2                               ! 2nd non-time dimension size (1 if absent)
      integer, intent(in) :: dimsiz3                               ! 3rd non-time dimension size (1 if absent)
      double precision, intent(in) :: arr(dimsiz1,dimsiz2,dimsiz3) ! array to be analyzed
      real(r8), intent(in) :: fillvalue                            ! fillvalue (if applicable)
      logical, intent(in) :: doprint                               ! whether or not to print requested data
      character(len=*), intent(in) :: dimnames(3)                  ! dimension names

      type(varinfo), intent(out) :: varstats                       ! stats will be put into this struct
!
! Local workspace
!
      integer :: i, j, k        ! spatial indices
      integer :: nval           ! number of array values
      integer :: imax           ! max value position indx 1st dimension
      integer :: jmax           ! max value position indx 2nd dimension
      integer :: kmax           ! max value position indx 3rd dimension
      integer :: imin           ! min value position indx 1st dimension
      integer :: jmin           ! min value position indx 2nd dimension
      integer :: kmin           ! min value position indx 3rd dimension
      integer :: nfill          ! number of fill values

      real(r8) :: avg           ! average field value
      real(r8) :: arrmax        ! max field value
      real(r8) :: arrmin        ! min field value

      logical :: chkfillvalue   ! whether to check for fillvalue
      logical :: validv         ! whether valid (i.e. non-fill) values exist
!
! Initialize stats
!
      arrmax = -huge (arrmax)
      arrmin = +huge (arrmax)
      avg    = 0.
      nval   = 0
      imax   = uninit_int
      jmax   = uninit_int
      kmax   = uninit_int
      imin   = uninit_int
      jmin   = uninit_int
      kmin   = uninit_int
      nfill  = 0
!
! Fillvalue will contain uninit_r8 if it does not apply to this variable
!
      chkfillvalue = fillvalue /= uninit_r8

      do k=1,dimsiz3
         do j=1,dimsiz2
            do i=1,dimsiz1
!
! First print data if so requested
!
               if (doprint) then
                  if (i >= iprs .and. i <= ipre .and. &
                      (dimsiz2 == 1 .or. (j >= jprs .and. j <= jpre)) .and. &
                      (dimsiz3 == 1 .or. (k >= kprs .and. k <= kpre))) then
                     write(6,100)i, j, k, arr(i,j,k)
100                  format(' i,j,k=',3i5,' arr=',1p,e23.15)
                  end if
               end if

               validv = .not. chkfillvalue .or. (chkfillvalue .and. arr(i,j,k) /= fillvalue)
               if (validv) then
                  avg = avg + abs (arr(i,j,k))
                  nval = nval + 1
                  if (arr(i,j,k) > arrmax) then
                     arrmax = arr(i,j,k)
                     imax   = i
                     jmax   = j
                     kmax   = k
                  end if
               
                  if (arr(i,j,k) < arrmin) then
                     arrmin = arr(i,j,k)
                     imin   = i
                     jmin   = j
                     kmin   = k
                  end if
               else
                  nfill = nfill + 1
               end if
            end do
         end do
      end do
!
! Copy output to varstats
!
      varstats%name         = name
      varstats%dimnames(:)  = dimnames(:)
      varstats%npossible    = dimsiz1*dimsiz2*dimsiz3
      varstats%nfill        = nfill
      if (nval > 0) then
         varstats%imax      = imax
         varstats%jmax      = jmax
         varstats%kmax      = kmax
         varstats%imin      = imin
         varstats%jmin      = jmin
         varstats%kmin      = kmin
         varstats%avg       = avg / nval
         varstats%arrmax    = arrmax
         varstats%arrmin    = arrmin
      else
         varstats%imax      = 0
         varstats%jmax      = 0
         varstats%kmax      = 0
         varstats%imin      = 0
         varstats%jmin      = 0
         varstats%kmin      = 0
         varstats%avg       = fillvalue
         varstats%arrmax    = fillvalue
         varstats%arrmin    = fillvalue
      end if

      return
   end subroutine gather_stats

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine gather_comparison_stats (name, ndims, dimsiz1, dimsiz2, dimsiz3, &
                                       arr, fillvalue, dimnames, diffstats, dodiffer)
!-------------------------------------------------------------------------------------------
! Purpose: Gather difference stats for the 2 instances of this variable. Stats include
!          number of diffs, largest absolute diff, largest relative diff, position 
!          info, and rms diff.
!
! Method:  Check only those points where _FillValue does not apply. Copy output into
!          diffstats. Also, set dodiffer to indicate whether any diffs were found
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      character(len=*), intent(in) :: name                           ! variable name
      integer, intent(in) :: ndims                                   ! number of dimensions                     
      integer, intent(in) :: dimsiz1                                 ! 1st non-time dimension size              
      integer, intent(in) :: dimsiz2                                 ! 2nd non-time dimension size (1 if absent)
      integer, intent(in) :: dimsiz3                                 ! 3rd non-time dimension size (1 if absent)
!
! Declare arrays as double precision not r8 for netcdf convenience
!
      double precision, intent(in) :: arr(dimsiz1,dimsiz2,dimsiz3,2) ! array to be analyzed
      real(r8), intent(in) :: fillvalue(2)                           ! fillvalue (if applicable)
      character(len=*), intent(in) :: dimnames(3)                    ! dimension names

      type(diffinfo), intent(out) :: diffstats                       ! difference stats will be put into this struct
      logical, intent(out) :: dodiffer                               ! flag indicates whether any differences exist
!
! Local workspace
!
      integer :: i, ii, j, k                   ! spatial indices
      integer :: isub(1)                       ! i-index
      integer :: ndiff                         ! number of diffs
      integer :: npossible                     ! number of field values
      integer :: indx(dimsiz1)                 ! indices where there are diffs
      integer :: nval                          ! intermediate number of diffs
      integer :: imax                          ! i-index of max diff
      integer :: jmax                          ! j-index of max diff
      integer :: kmax                          ! k-index of max diff
      integer :: irmax                         ! i-index of max relative diff
      integer :: jrmax                         ! j-index of max relative diff
      integer :: krmax                         ! k-index of max relative diff

      real(r8) :: meanvertwt                   ! mean vertwt between 2 cases
      real(r8) :: diff(dimsiz1)                ! difference
      real(r8) :: dmax                         ! max difference
      real(r8) :: rdiff(dimsiz1)               ! relative difference
      real(r8) :: rdmax                        ! max relative difference
      real(r8) :: dvals(2)                     ! array values where max diff occurs
      real(r8) :: rdvals(2)                    ! array values where max relative diff occurs
      real(r8) :: rms                          ! RMS difference
      real(r8) :: denom                        ! denominator of expression
      real(r8) :: rdbar                        ! used in computing mean relative differenc
      real(r8) :: rdlnbar                      ! used in computing mean number of digits that match

      logical :: chkfillvalue1, chkfillvalue2  ! whether to check for fillvalue
      logical :: validv1, validv2              ! whether an element equals fillvalue or not

      ndiff     = 0
      npossible = 0
      imax      = uninit_int
      jmax      = uninit_int
      kmax      = uninit_int
      irmax     = uninit_int
      jrmax     = uninit_int
      krmax     = uninit_int
      dvals(:)  = uninit_r8
      rdvals(:) = uninit_r8
      dmax      = -huge (dmax)
      rdmax     = -huge (dmax)
      rdbar     = 0.
      rdlnbar   = 0.
      rms       = 0.

      chkfillvalue1 = fillvalue(1) /= uninit_r8
      chkfillvalue2 = fillvalue(2) /= uninit_r8

      do k=1,dimsiz3
         do j=1,dimsiz2
            nval = 0
            do i=1,dimsiz1
!
! First print data if so requested
!
               if (i >= iprs .and. i <= ipre .and. j >= jprs .and. j <= jpre .and. &
                   (dimsiz3 == 1 .or. k >= kprs .and. k <= kpre)) then
                  write(6,100)i, j, k, arr(i,j,k,1), arr(i,j,k,2), arr(i,j,k,1) - arr(i,j,k,2)
100               format(' i,j,k=',3i5,' arr1,arr2,diff=',1p,3e23.15)
               end if

               npossible = npossible + 1
               validv1 = .not. chkfillvalue1 .or. (chkfillvalue1 .and. arr(i,j,k,1) /= fillvalue(1))
               validv2 = .not. chkfillvalue2 .or. (chkfillvalue2 .and. arr(i,j,k,2) /= fillvalue(2))

               if (validv1 .and. validv2) then
                  diff(i) = abs (arr(i,j,k,1) - arr(i,j,k,2))
                  rms = rms + diff(i)**2
                  if (diff(i) /= 0.) then
                     nval = nval + 1
                     indx(nval) = i
                  end if
               else
!
! Setting diff to zero here means locations with fillvalue in one file and valid data
! in the other file will be treated as not different.  This is probably not the best
! way to handle things.
!
                  diff(i) = 0.
               end if
            end do
!
! Save max diff info
! 
            if (nval > 0) then                     ! at least 1 difference
               ndiff = ndiff + nval
               isub = maxloc (diff(:))
               if (diff(isub(1)) > dmax) then      ! Save max diff info
                  dmax = diff(isub(1))
                  dvals(1) = arr(isub(1),j,k,1)
                  dvals(2) = arr(isub(1),j,k,2)
                  imax = isub(1)
                  jmax = j
                  kmax = k
               end if
!
! Save max relative diff info
!
               rdiff(:) = 0.
               do ii=1,nval                              ! Compute relative diffs
                  i = indx(ii)
                  denom = max (abs (arr(i,j,k,1)), abs (arr(i,j,k,2)))
                  rdiff(i) = diff(i)/(2.*denom)
                  rdbar = rdbar + rdiff(i)
                  rdlnbar = rdlnbar - log10 (rdiff(i))
               end do
               isub = maxloc (rdiff(:))
            
               if (rdiff(isub(1)) > rdmax) then          ! Save max relative diff info
                  rdmax = rdiff(isub(1))
                  rdvals(1) = arr(isub(1),j,k,1)         ! Save values & indices
                  rdvals(2) = arr(isub(1),j,k,2)
                  irmax = isub(1)
                  jrmax = j
                  krmax = k
               end if
            end if          ! nval > 0
         end do             ! j=1,dimsiz2
      end do                ! k=1,dimsiz3
!
! Final normalization
!
      rms   = sqrt (rms/npossible)
      rdbar = rdbar/npossible
!
! Copy values into struct
!
      diffstats%name        = name
      diffstats%dimnames(:) = dimnames(:)
      diffstats%dmax        = dmax
      diffstats%rdmax       = rdmax
      diffstats%npossible   = npossible
      diffstats%ndiff       = ndiff
      diffstats%imax        = imax
      diffstats%jmax        = jmax
      diffstats%kmax        = kmax
      diffstats%irmax       = irmax
      diffstats%jrmax       = jrmax
      diffstats%krmax       = krmax
      diffstats%dvals(:)    = dvals(:)
      diffstats%rdvals(:)   = rdvals(:)
      diffstats%rms         = rms
      diffstats%rdbar       = rdbar
      diffstats%rdlnbar     = rdlnbar

      dodiffer = ndiff > 0

      return
   end subroutine gather_comparison_stats

end module gatherstats
