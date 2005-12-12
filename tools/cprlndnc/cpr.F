subroutine cpr (nlon, nlat, nlev, numcases, lnorm, verbose)
  !
  ! $Id$
  !
  use precision
  use header
  use stats
  use nldat

  implicit none

  include 'netcdf.inc'
  !
  ! Input arguments
  !
  integer, intent(in) :: nlon, nlat, nlev, numcases

  logical lnorm            ! Whether to print lvl-by-lvl L2 & L-INF norm stats
  logical verbose
  !
  ! Local workspace (static or stack)
  !
  real(r8), parameter :: missing_value = 1.e+36    ! hard-wire for now

  character*(nf_max_name) name

  integer nvars            ! number of variables on tape
  integer natts
  integer n,n2             ! variable index
  integer ndims(2)
  integer idim
  integer xtype(2)
  integer dimids(nf_max_dims,2)
  integer i, j, k          ! spatial indices
  integer numlev,numlev2   ! number of levels for each field
  integer isub, nval, ii
  integer numvalid(2)
  integer itime, itime2    ! loop index over time
  integer ret              ! return code
  integer, dimension(3) :: start2d, count2d

  real(r8) wsum
  real(r8) pdelbar
  real(r8) denom
  real(r8) rd
  real(r8), parameter :: timeepsilon = 1.e-9
  real(r8) :: timediff

  logical twocases         ! True => 2 h-tapes are being analyzed
  logical found
  !
  ! Local workspace (dynamic)
  !
  real(r8), dimension(nlon,nlat,nlev) :: arr, arr2
  real(r8) diff(nlon)             ! difference field
  real(r8) pmid(nlon,nlev)
  real(r8) rdiff(nlon)            ! relative difference field

  integer indx(nlon)          ! indices of non-zero diffs
  !
  ! Externals
  !
  integer ismax, fillarr, lenchr
  external ismax, fillarr, lenchr
  !
  !------------------------------------------------------------------------------
  !
  ! Define start, count arrays
  !
  start2d(:) = 1
  count2d(1) = nlon
  count2d(2) = nlat
  count2d(3) = 1

  twocases = numcases.eq.2
  !
  ! Find number of variables on the 1st tape
  !
  call wrap_inq_nvars (ncid(1), nvars)
  !
  ! Big loop over number of time slices.  Not written as explicit "do" since
  ! may need to adjust time index in the middle of the loop to match time
  ! levels on tapes.
  !
  itime = 1
  itime2 = 1
10 continue
  !
  ! Check equality of time steps
  ! For lnd model just skip since time variable not available
  !
  if (.false.) then
     if (matchts .and. twocases) then
        timediff = abs (time(itime) - time2(itime2))
        do while (timediff .gt. timeepsilon)
           write(6,'(1x,a,f12.4,/1x,a,f12.4)')&
                'Current times differ:time1=', time(itime), &
                '                     time2=', time2(itime2)
           if (time(itime).lt.time2(itime2)) then
              itime = itime + 1
              if (itime.gt.ntime(1)) then
                 write(6,*)'End of tape 1 reached: stopping'
                 stop 0
              end if
           else
              itime2 = itime2 + 1
              if (itime2.gt.ntime(2)) then
                 write(6,*)'End of tape 2 reached: stopping'
                 stop 0
              end if
           end if
           timediff = abs (time(itime) - time2(itime2))
        end do
     end if
  end if

  if (twocases) then
     write(6,800)case(1)(1:lenchr(case(1))), case(2)(1:lenchr(case(2))),  &
          title(1)(1:lenchr(title(1))), title(2)(1:lenchr(title(2))), &
          date_written(itime,1), date_written(itime2,2), &
          time_written(itime,1), time_written(itime2,2), &
          itime, itime2, &
          ntime(1), ntime(2), &
          nsteph(itime,1), nsteph(itime2,2), &
          ncdate(itime,1), ncdate(itime2,2), &
          ncsec(itime,1),  ncsec(itime2,2)

  else                               ! One tape only.

     write(6,801)case(1)(1:lenchr(case(1))),  &
          title(1)(1:lenchr(title(1))), &
          date_written(itime,1), &
          time_written(itime,1), &
          itime, &
          ntime(1), &
          nsteph(itime,1), &
          ncdate(itime,1), &
          ncsec(itime,1)
  end if                    ! twocases
  !
  ! Read in and check the data fields sequentially
  !
  if (iprs.ne.0) then
     write(6,'(/,a,i4,a,i4, a,i3,a,i3, a,i3,a,i3)') &
          ' SUMMARY OF FIELD INFORMATION FOR LONGITUDE INDICES ', &
          iprs,' through ',ipre,', LATITUDE INDICES ', &
          latprs,' through ',latpre,', LEVELS ',kprs,' through ',kpre

     write(6,'(a,t10,a,t14,a,t18,a,t23,a,t35,a,t47,a,t70,a)') &
          ' FIELD','LON','LAT','LEV','DIFF','RDIFF','CASE1','CASE2'
  end if

  wsum = 0.
  if (twocases) then
     do j=1,nlat
        do i= 1,nlon
           wsum = wsum + area(i,j,1)
        end do
     end do
  end if
  !
  ! Loop over all variables on the 1st tape
  !
  do n=1,nvars

     dimids(:,1) = -1
     call wrap_inq_var (ncid(1), n, name, xtype(1), ndims(1), dimids(1,1), natts)

     if (verbose) then
        write(6,*) 'processing ', name(1:lenchr(name))
     end if
     !
     ! Determine number of levels for this variable
     !
     found = .false.
     numlev = 1
     do idim=1,ndims(1)
        if (dimids(idim,1).eq.levdimid(1)) then
           numlev = nlev
           if (found) then
              write(6,*) name(1:lenchr(name)),' has 2 level dims: skipping'
              cycle
           end if
           found = .true.
        end if
     end do

     ret = fillarr (ncid(1), n, dimids(1,1), londimid(1), latdimid(1), &
          levdimid(1), nlon, nlat, &
          numlev, arr, itime)
     !
     ! Skip to the next field if fillarr failed
     !
     if (ret.lt.0) cycle
     !
     ! Initialize statistics, then modify for this field
     !
     call initstats

     xmx(1) = -1.d99
     xmn(1) = +1.d99
     numvalid(:) = 0
     do k=1,numlev
        do j=1,nlat
           do i=1,nlon
              if (arr(i,j,k) < missing_value) then
                 numvalid(1) = numvalid(1) + 1
                 xbar(1) = xbar(1) + abs(arr(i,j,k))
                 if (arr(i,j,k).gt.xmx(1)) then
                    xmx(1) = arr(i,j,k)
                    imx(1) = i
                    jmx(1) = j
                    kmx(1) = k
#ifdef DEBUG
                    if (xmx(1) > 1.e35) then
                       write(6,*)'i,j,k,xmx=',i,j,k,xmx(1)
                       stop
                    end if
#endif
                 end if

                 if (arr(i,j,k).lt.xmn(1)) then
                    xmn(1) = arr(i,j,k)
                    imn(1) = i
                    jmn(1) = j
                    kmn(1) = k
                 end if
              end if
           end do
        end do
     end do

     if (twocases) then
        ret = nf_inq_varid (ncid(2), name, n2)
        if (ret.ne.NF_NOERR) then
           if (verbose) then
              write(6,*) name, ' not found on tape 2'
           end if
           cycle
        end if

        dimids(:,2) = -1
        call wrap_inq_var (ncid(2), n2, name, xtype(2), ndims(2), &
             dimids(1,2), natts)

        if (xtype(1).ne.xtype(2)) then
           write(6,*)'NOTE: type variables do not match for field ', &
                name(1:lenchr(name))
        end if

        if (xtype(1).eq.NF_FLOAT .or. xtype(2).eq.NF_FLOAT) then
           write(6,*)'NOTE: At least 1 tape has 32-bit values for field ', &
                name(1:lenchr(name))
        end if
        !
        ! Determine number of levels for this variable
        !
        found = .false.
        numlev2 = 1
        do idim=1,ndims(2)
           if (dimids(idim,2).eq.levdimid(2)) then
              numlev2 = nlev
              if (found) then
                 write(6,*) name(1:lenchr(name)),' has 2 level dims: skipping'
                 cycle
              end if
              found = .true.
           end if
        end do

        if (numlev2.ne.numlev) then
           write(6,*) name(1:lenchr(name)), ' numlev mismatch: skipping'
           cycle
        end if

        ret = fillarr (ncid(2), n2, dimids(1,2), londimid(2), latdimid(2), &
             levdimid(2), nlon, nlat, &
             numlev, arr2, itime)
        !
        ! Skip to the next field if fillarr failed
        !
        if (ret.lt.0) then
           write(6,*)'Tape 2 does not have field ',name,': skipping'
           cycle
        end if

        xmx(2) = -1.d99
        xmn(2) = +1.d99
        do k=1,numlev
           do j=1,nlat

              ! No need to mask missing_value since difference will be zero

              diff(:) = abs (arr(:,j,k) - arr2(:,j,k))

              do i=1,nlon
                 diffmw = diffmw + (diff(i)**2)*area(i,j,1)
              end do

              do i=1,nlon
                 if (arr(i,j,k) < missing_value .and. arr2(i,j,k) < missing_value) then
                    npos = npos + 1
                    rms = rms + diff(i)**2
                    xbar(2) = xbar(2) + abs(arr2(i,j,k))
                    if (arr2(i,j,k).gt.xmx(2)) then
                       xmx(2) = arr2(i,j,k)
                       imx(2) = i
                       jmx(2) = j
                       kmx(2) = k
                    end if
#ifdef DEBUG
                    if (xmx(2) > 1.e35) then
                       write(6,*)'i,j,k,xmx=',i,j,k,xmx(2)
                       write(6,*)'arr,arr2=',arr(i,j,k),arr2(i,j,k)
                       stop
                    end if
#endif
                    if (arr2(i,j,k).lt.xmn(2)) then
                       xmn(2) = arr2(i,j,k)
                       imn(2) = i
                       jmn(2) = j
                       kmn(2) = k
                    end if
                 end if
              end do

              if (iprs.ne.0) then
                 if (numlev.eq.1 .or. (k.ge.kprs .and. k.le.kpre)) then
                    if (j.ge.latprs .and. j.le.latpre) then
                       do i=iprs,ipre
                          denom = 0.5*(abs(arr(i,j,k)) + abs(arr2(i,j,k)))
                          rd = 0.
                          if (denom.ne.0.) rd = diff(i)/denom
                          write(6,910) name,i,j,k,diff(i),rd,arr(i,j,k),arr2(i,j,k)
                       end do
                    end if
                 end if
              end if
              !
              ! Test on half of difference field rather than full field due to 0.5 factor
              ! used in computation of "denom" later.
              ! 
              nval = 0
              do i=1,nlon
                 if (diff(i).ne.0.) then
                    nval = nval + 1
                    indx(nval) = i
                 end if
              end do

              if (nval.gt.0) then                 ! at least 1 difference
                 ndif = ndif + nval
                 isub = ismax(nlon, diff, 1)
                 if (diff(isub).gt.difmx) then    ! Save max diff info
                    difmx = diff(isub)
                    dmxsv(1) = arr(isub,j,k)      ! Save values
                    dmxsv(2) = arr2(isub,j,k)
                    idmxsv = isub
                    jdmxsv = j
                    kdmxsv = k
                 end if

                 rdiff(:) = 0.
                 do ii=1,nval                     ! Compute relative diffs
                    i = indx(ii)
                    denom = max(abs(arr(i,j,k)), abs(arr2(i,j,k)))
                    rdiff(i) = diff(i)/(2.*denom)
                    rdbar = rdbar + rdiff(i)
                    rdlnbar = rdlnbar - log10(rdiff(i))
                 end do
                 isub = ismax(nlon,rdiff,1)

                 if (rdiff(isub).gt.rdifmx) then  ! Save max relative diff info
                    rdifmx = rdiff(isub)
                    rdmxsv(1) = arr(isub,j,k)     ! Save values & indices
                    rdmxsv(2) = arr2(isub,j,k)
                    irdmxsv = isub
                    jrdmxsv = j
                    krdmxsv = k
                 end if
              end if

           end do
        end do
     end if
     !
     ! Print stats
     !
     call printstats (twocases, name, numlev, nlev, &
                      nlon, nlat, wsum)
  end do  ! n=1,nvars
  !
  ! Proceed to next time slice if more are available
  !
  itime = itime + 1
  if (itime.gt.ntime(1)) then
     write(6,*)'End of tape 1 reached: stopping'
     stop 0
  end if

  if (twocases) then
     itime2 = itime2 + 1
     if (itime2.gt.ntime(2)) then
        write(6,*)'End of tape 2 reached: stopping'
        stop 0
     end if
  end if

  goto 10

800 format(' CASE:',2(A,1X),/, &
       ' CASE1 TITLE:',a80,/, ' CASE2 TITLE:',a80,/, &
       ' date_written: ',2(A8,1X),/, &
       ' time_written: ',2(A8,1X),/, &
       ' MFILH:  ',2I9,/,' MFILTH: ',2I9,/,' NSTEPH: ',2I9,/, &
       ' NCDATE: ',2I9,/, &
       ' NCSEC:  ',2I9/)

801 format(' CASE1:',a,/, &
       ' CASE1 TITLE:',a80,/, &
       ' date_written: ',1(A8,1X),/, &
       ' time_written: ',1(A8,1X),/, &
       ' MFILH:  ',1I9,/,' MFILTH: ',1I9,/,' NSTEPH: ',1I9,/, &
       ' NCDATE: ',1I9,/, &
       ' NCSEC:  ',1I9,/)

910 format(1x,a8,'(',i3,',',i2,',',i2,')=',1p,2e12.4,2e23.15)

end subroutine cpr

integer function fillarr (ncid, varid, dimids, londimid, latdimid, &
                          levdimid, nlon, nlat, &
                          numlev, arr, itime)
  use precision

  implicit none

  include 'netcdf.inc'

  integer :: ncid, varid
  integer :: dimids(NF_MAX_DIMS)
  integer :: londimid, latdimid, levdimid
  integer :: nlon, nlat, numlev
  integer :: itime

  real(r8) :: arr(nlon,nlat,numlev)
!
! Local workspace
!
  integer j, k
  integer start(4), count(4)

  real(r8), allocatable :: arrxzy(:,:,:)

  fillarr = 0       ! Initialize return value to success

  if (dimids(1).eq.londimid .and. &
      dimids(2).eq.latdimid .and. &
      (dimids(3).eq.levdimid)) then

    start(1) = 1
    start(2) = 1
    start(3) = 1
    start(4) = itime

    count(1) = nlon
    count(2) = nlat
    count(3) = numlev
    count(4) = 1

    call wrap_get_vara_realx (ncid, varid, start, count, arr)

  else if (dimids(1).eq.londimid .and. &
           (dimids(2).eq.levdimid) .and. &
           dimids(3).eq.latdimid ) then

    start(1) = 1
    start(2) = 1
    start(3) = 1
    start(4) = itime

    count(1) = nlon
    count(2) = numlev
    count(3) = nlat
    count(4) = 1

    allocate(arrxzy(nlon,numlev,nlat))
    call wrap_get_vara_realx (ncid, varid, start, count, arrxzy)
    do k=1,numlev
      do j=1,nlat
        arr(:,j,k) = arrxzy(:,k,j)
      end do
    end do
    deallocate(arrxzy)

  else if (dimids(1).eq.londimid .and. dimids(2).eq.latdimid) then

    start(1) = 1
    start(2) = 1
    start(3) = itime
    start(4) = 0

    count(1) = nlon
    count(2) = nlat
    count(3) = 1
    count(4) = 0

    call wrap_get_vara_realx (ncid, varid, start, count, arr)

  else

    fillarr = -1

  end if

  return
end function fillarr

#if ( ! defined TIMING )
subroutine t_startf (xxx)
  character*(*) xxx
  return
end subroutine t_startf

subroutine t_stopf (xxx)
  character*(*) xxx
  return
end subroutine t_stopf

subroutine t_prf (xxx)
  integer xxx
  return
end subroutine t_prf
#endif
