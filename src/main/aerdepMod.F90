module aerdepMOD

!-----------------------------------------------------------------------
! This entire module will be removed....................
!BOP
!
! !MODULE: aerdepMod
!
! !DESCRIPTION:
! read an interpolate aerosol deposition data
!
! !USES:
  use shr_kind_mod,    only : r8 => shr_kind_r8
  use abortutils,      only : endrun
  use clm_varctl,      only : scmlat,scmlon,single_column
  use clm_varctl,      only : iulog
  use clm_varcon     , only : secspday
  use perf_mod,        only : t_startf, t_stopf
  use clm_varctl,      only : set_caerdep_from_file, set_dustdep_from_file
  use spmdMod,         only : masterproc, mpicom, MPI_REAL8, MPI_INTEGER
  use fileutils      , only : getfil
  use ncdio          , only : check_ret,ncd_iolocal
  use shr_sys_mod    , only : shr_sys_flush


!
! !PUBLIC TYPES:
  implicit none

  private

! !INCLUDES:
  include 'netcdf.inc'
  
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: interpMonthlyAerdep   ! interpolate monthly aerosol deposition data
  public :: aerdepini             ! aerosol deposition initialization

!
! !REVISION HISTORY:
! Created by Mark Flanner, 
!   based on vegetation interpolation schemes in STATICEcosystemDynMod
!    2009-Apr-17 B. Kauffman -- added multi-year time series functionality
!
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: readMonthlyAerdep       ! read monthly aerosol deposition data for two months

!
! !PRIVATE TYPES:
!EOP

  real(r8), save, private, allocatable :: bcphiwet2t(:,:)
  real(r8), save, private, allocatable :: bcphidry2t(:,:)
  real(r8), save, private, allocatable :: bcphodry2t(:,:)
  real(r8), save, private, allocatable :: ocphiwet2t(:,:)
  real(r8), save, private, allocatable :: ocphidry2t(:,:)
  real(r8), save, private, allocatable :: ocphodry2t(:,:)
  real(r8), save, private, allocatable :: dstx01wd2t(:,:)
  real(r8), save, private, allocatable :: dstx01dd2t(:,:)
  real(r8), save, private, allocatable :: dstx02wd2t(:,:)
  real(r8), save, private, allocatable :: dstx02dd2t(:,:)
  real(r8), save, private, allocatable :: dstx03wd2t(:,:)
  real(r8), save, private, allocatable :: dstx03dd2t(:,:)
  real(r8), save, private, allocatable :: dstx04wd2t(:,:)
  real(r8), save, private, allocatable :: dstx04dd2t(:,:)

  integer             ,save :: nt      ! size of time(:) array
  real(r8),allocatable,save :: time(:) ! data time, elapsed days since 0000-01-01 0s
  real(r8),parameter  :: daysPerYear = 365.0_r8 ! days per year

  integer,parameter :: debug = 1 ! internal debug level

!================================================================================

contains

!================================================================================
!BOP
!
! !IROUTINE: aerdepini
!
! !INTERFACE:
  subroutine aerdepini()
!
! !DESCRIPTION:
! Dynamically allocate memory and set to signaling NaN.
!
! !USES:
    use nanMod         , only : nan
    use decompMod      , only : get_proc_bounds
    use clm_varctl     , only : faerdep
    use shr_ncread_mod , only : shr_ncread_tCoord
    use shr_cal_mod    , only : shr_cal_date2eday
!
! !ARGUMENTS:
    implicit none

!
! !REVISION HISTORY:
!    2009-Apr-17 B. Kauffman -- added multi-year time series functionality
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: ier         ! error code
    integer :: begg,endg   ! local beg and end p index

    character(256) :: locfn          ! local file name
    integer :: ncid,dimid,varid      ! input netCDF id's

    integer,allocatable :: cdate(:)  ! calendar date yyyymmdd
    integer,allocatable :: eday(:)   ! elapsed days since 0000-01-01
    integer,allocatable :: secs(:)   ! elapsed secs within current date
    integer             :: n         ! loop index
    integer             :: m1,m2     ! month 1, month 2

    integer, parameter :: ndaypm(12) = &
         (/31,28,31,30,31,30,31,31,30,31,30,31/) !days per month

    character(*),parameter :: subName =   '(aerdepini) '
    character(*),parameter :: F00     = "('(aerdepini) ',4a)"
    character(*),parameter :: F01     = "('(aerdepini) ',a,4f13.3)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

    call get_proc_bounds(begg=begg,endg=endg)

    if ( set_caerdep_from_file )then
       ier = 0
       if(.not.allocated(bcphiwet2t)) then
          allocate(bcphiwet2t(begg:endg,2))
          allocate(bcphidry2t(begg:endg,2))
          allocate(bcphodry2t(begg:endg,2))
          allocate(ocphiwet2t(begg:endg,2))
          allocate(ocphidry2t(begg:endg,2))
          allocate(ocphodry2t(begg:endg,2))
       endif
          
       if (ier /= 0) call endrun( 'aerdepini allocation error' )

        bcphiwet2t(begg:endg,1:2) = nan
        bcphidry2t(begg:endg,1:2) = nan
        bcphodry2t(begg:endg,1:2) = nan
        ocphiwet2t(begg:endg,1:2) = nan
        ocphidry2t(begg:endg,1:2) = nan
        ocphodry2t(begg:endg,1:2) = nan
    end if

    if ( set_dustdep_from_file )then

       if(.not.allocated(dstx01wd2t)) then
          allocate(dstx01wd2t(begg:endg,2))
          allocate(dstx01dd2t(begg:endg,2))
          allocate(dstx02wd2t(begg:endg,2))
          allocate(dstx02dd2t(begg:endg,2))
          allocate(dstx03wd2t(begg:endg,2))
          allocate(dstx03dd2t(begg:endg,2))
          allocate(dstx04wd2t(begg:endg,2))
          allocate(dstx04dd2t(begg:endg,2))
       endif

       if (ier /= 0) call endrun( 'aerdepini allocation error' )

       dstx01wd2t(begg:endg,1:2) = nan
       dstx01dd2t(begg:endg,1:2) = nan
       dstx02wd2t(begg:endg,1:2) = nan
       dstx02dd2t(begg:endg,1:2) = nan
       dstx03wd2t(begg:endg,1:2) = nan
       dstx03dd2t(begg:endg,1:2) = nan
       dstx04wd2t(begg:endg,1:2) = nan
       dstx04dd2t(begg:endg,1:2) = nan

    end if

   !----------------------------------------------------------------------------
   ! read time axis from data file
   !----------------------------------------------------------------------------

   if ( masterproc ) write(iulog,F00) "Read time axis from ",trim(faerdep)

   if (masterproc) then
      call getfil(faerdep, locfn, 0)
      call check_ret(nf_open      (locfn, 0, ncid     ), subname)
      call check_ret(nf_inq_dimid (ncid, 'time', dimid), subname)
      call check_ret(nf_inq_dimlen(ncid, dimid, nt    ), subname)
   endif 
   call mpi_bcast (nt, 1, MPI_INTEGER, 0, mpicom, ier)
   if ( masterproc ) write(iulog,*)  subName,"number of time samples in file = ",nt
   allocate(time(nt))

   if ( nt < 12) then 
      write(iulog,F00)     "ERROR: input data must have at least 12 months of data"
      call endrun(subName//"ERROR: input data must have at least 12 months of data")
   else if ( nt == 12) then 
      !--- old format => create time axis here ---
      if ( masterproc ) write(iulog,F00)  'nt = 12: assume old format without a CF-1.0 time axis'
      n = 365
      time( 1) = n + float(ndaypm( 1))/2.0_r8 ; n = n + ndaypm( 1)
      time( 2) = n + float(ndaypm( 2))/2.0_r8 ; n = n + ndaypm( 2)
      time( 3) = n + float(ndaypm( 3))/2.0_r8 ; n = n + ndaypm( 3)
      time( 4) = n + float(ndaypm( 4))/2.0_r8 ; n = n + ndaypm( 4)
      time( 5) = n + float(ndaypm( 5))/2.0_r8 ; n = n + ndaypm( 5)
      time( 6) = n + float(ndaypm( 6))/2.0_r8 ; n = n + ndaypm( 6)
      time( 7) = n + float(ndaypm( 7))/2.0_r8 ; n = n + ndaypm( 7)
      time( 8) = n + float(ndaypm( 8))/2.0_r8 ; n = n + ndaypm( 8)
      time( 9) = n + float(ndaypm( 9))/2.0_r8 ; n = n + ndaypm( 9)
      time(10) = n + float(ndaypm(10))/2.0_r8 ; n = n + ndaypm(10)
      time(11) = n + float(ndaypm(11))/2.0_r8 ; n = n + ndaypm(11)
      time(12) = n + float(ndaypm(12))/2.0_r8

      if ( masterproc ) write(iulog,F01) 'first time series cal date = ',10100.0+ndaypm( 1)/2.0
      if ( masterproc ) write(iulog,F01) 'last  time series cal date = ',11200.0+ndaypm(12)/2.0
   else 
      !--- new format => get time axis from data file ---
      if (masterproc) then
         write(iulog,F00)  "nt > 12: assume new format with a CF-1.0 time axis"
         allocate(eday(nt),cdate(nt),secs(nt))

         call shr_ncread_tCoord(locfn,"time",cdate,secs,ier)

         !--- convert date & secs to elapsed days since 0000-01-01 0sec ---
         do n = 1,nt
            call shr_cal_date2eday(cdate(n),eday(n))
            time(n) = float(eday(n)) + float(secs(n))/secspday
            if (debug>2) write(iulog,*) subName,"n,cdate,secs,time",n,cdate(n),secs(n),time(n)
         end do
         write(iulog,F01) 'first time series cal date = ',cdate( 1) + secs( 1)/secspday
         write(iulog,F01) 'last  time series cal date = ',cdate(nt) + secs(nt)/secspday

         !--- requirement: 1 time sample per month & don't skip months ---
         do n = 2,nt
            m1 = mod( cdate(n-1)/100, 100) + 1 ! get month, add one
            m2 = mod( cdate(n  )/100, 100) 
            if (m1==13) m1 = 1
            if (m1 .ne. m2) then
                write(iulog,*  )     subname,"cdate(n-1),cdate(n)=",cdate(n-1),cdate(n)
                write(iulog,*  )     subname,"m1        ,m2      =",m1,m2
                write(iulog,F00)     "ERROR: data must have exactly one sample per month"
                call endrun(subName//"ERROR: data must have exactly one sample per month")
            end if
            if (m1 == 1) then
                m1 = cdate(n-1)/10000 + 1 ! get year, add one
                m2 = cdate(n  )/10000 
                if (m1 .ne. m2) then
                   write(iulog,*  )     subname,"cdate(n-1),cdate(n)=",cdate(n-1),cdate(n)
                   write(iulog,*  )     subname,"m1        ,m2      =",m1,m2
                   write(iulog,F00)     "ERROR: input data skipped a year"
                   call endrun(subName//"ERROR: input data skipped a year")
               end if
            end if
         end do

         !--- requirement: 1st month is Jan, last is Dec ---
         m1 = mod( cdate(1 )/100, 100) ! get month, add one
         m2 = mod( cdate(nt)/100, 100) 
         if (m1.ne.1  .or.  m2.ne.12) then
             write(iulog,F00)     "ERROR: first month must be Jan, last Dec"
             write(iulog,*  )     subName,"ERROR: first & last month = ",m1,m2
             call endrun(subName//"ERROR: first month must be Jan, last Dec")
         end if

         deallocate(eday,cdate,secs)
      endif 
      call mpi_bcast(time,size(time),MPI_REAL8,0,mpicom,ier)
   end if
   if ( masterproc )then
      write(iulog,F01) 'first time axis elapsed days since 0000-01-01 = ',time( 1)
      write(iulog,F01) 'last  time axis elapsed days since 0000-01-01 = ',time(nt)
      call shr_sys_flush(iulog)
   end if


  end subroutine aerdepini

!================================================================================
!BOP
!
! !IROUTINE: interpMonthlyAerdep
!
! !INTERFACE:
  subroutine interpMonthlyAerdep ()
!
! !DESCRIPTION:
! Determine if 2 new months of data are to be read.
!
! !USES:
    use clmtype
    use clm_atmlnd        , only : clm_a2l
    use clm_varctl        , only : faerdep
    use clm_time_manager  , only : get_curr_date, get_step_size, get_perp_date, is_perpetual
    use decompMod         , only : get_proc_bounds
    use shr_cal_mod       , only : shr_cal_ymd2eday
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!    2009-Apr-17 B. Kauffman -- added multi-year time series functionality
!  Adapted by Mark Flanner
!
!

!
! local pointers to implicit out arguments
!
    real(r8), pointer :: forc_aer(:,:)   ! aerosol deposition rate (kg/m2/s)

! !LOCAL VARIABLES:
!EOP
    real(r8):: timwt_aer(2)  ! time weights for month 1 and month 2 (aerosol deposition)
    integer :: kyr         ! year (0, ...) for nstep+1
    integer :: kmo         ! month (1, ..., 12)
    integer :: kda         ! day of month (1, ..., 31)
    integer :: ksec        ! seconds into current date for nstep+1
    real(r8):: dtime       ! land model time step (sec)
    integer :: g
    integer :: begg,endg                  ! beg and end local g index

    integer       :: n         ! counter to prevent infinite LB/UB search
    integer       :: edays     ! elapsed days since 0000-01-01 0s (excluding partial days)
    real(r8)      :: t         ! model time, elapsed days since 0000-01-01
    integer ,save :: nLB = 0   ! tLB = time(nLB)
    integer ,save :: nUB = 1   ! tUB = time(nUB)
    real(r8),save :: tLB =-1.0 ! upper bound time sample, model time is in [tLB,tUB]
    real(r8),save :: tUB =-2.0 ! lower bound time sample, model time is in [tLB,tUB]
    real(r8)      :: fUB,fLB   ! t-interp fracs for UB,LB
    logical ,save :: firstCallA = .true.   ! id 1st occurance of case A
    logical ,save :: firstCallB = .true.   ! id 1st occurance of case B
    logical ,save :: firstCallC = .true.   ! id 1st occurance of case C
    character(1)  :: case                  ! flags case A, B, or C
    logical       ::  readNewData          ! T <=> read new LB,UB data
    
    character(*),parameter :: subName =  '(interpMonthlyAerdep) '
    character(*),parameter :: F00    = "('(interpMonthlyAerdep) ',4a)"
    character(*),parameter :: F01    = "('(interpMonthlyAerdep) ',a,i4.4,2('-',i2.2),3f11.2,2i6,2x,2f6.3)"
    character(*),parameter :: F02    = "('(interpMonthlyAerdep) ',a,i4.4,2('-',i2.2),i7,'s ',f12.3)"

!-------------------------------------------------------------------------------
! WARNING: this is (and has always been) hard-coded to assume 365 days per year
!-------------------------------------------------------------------------------

    ! Determine necessary indices
    call get_proc_bounds(begg=begg,endg=endg)
    

    ! Assign local pointers to derived subtypes components (gridcell level)
    forc_aer     => clm_a2l%forc_aer
    
   !----------------------------------------------------------------------------
   ! get model time, convert units to days elapsed since 0000-01-01
   !----------------------------------------------------------------------------
    dtime = get_step_size()

    if ( is_perpetual() ) then
       call get_perp_date(kyr, kmo, kda, ksec, offset=int(dtime))
    else
       call get_curr_date(kyr, kmo, kda, ksec, offset=int(dtime))
    end if

   call shr_cal_ymd2eday(kyr,kmo,kda,edays) ! convert to elapsed days
   t = float(edays) + ksec/secspday
   if (masterproc .and. debug > 2) then
      write(iulog,F02) "model date, elapsed days = ", kyr,kmo,kda,ksec,t
   endif

   !----------------------------------------------------------------------------
   ! find input data LB & UB, time units are elapsed days since 0000-01-01
   !----------------------------------------------------------------------------

   CASE = "B"                    ! => interpolate within input time series
   if (t < time( 1) ) CASE = "A" ! => loop over 1st  year of input data
   if (t > time(nt) ) CASE = "C" ! => loop over last year of input data

   if ( case == "A" ) then 
      !--- CASE A: loop over first year of data ----------------------
      if ( firstCallA ) then
         if ( masterproc ) write(iulog,F00) "CASE A: loop over first year of data"
         if ( masterproc ) write(iulog,F01) "CASE A: model date, t, time(1) = ",kyr,kmo,kda,t,time(1)
         nLB = 0 ; tLB = -2.0
         nUB = 1 ; tUB = -1.0 ! forces search for new LB,UB
         firstCallA = .false.
      end if
      t = mod(t,daysPerYear) + daysPerYear ! CASE A: put t in year 1
      n = 0
      readNewData = .false.
      do while (t < tLB  .or.  tUB < t)
         readNewData = .true.
         !--- move tUB,tLB forward in time ---
         nLB = nLB + 1 ; if (nLB > 12) nLB = 1
         nUB = nLB + 1 ; if (nUB > 12) nUB = 1
         tLB = mod(time(nLB),daysPerYear) + daysPerYear ! set year to 1
         tUB = mod(time(nUB),daysPerYear) + daysPerYear 
         !--- deal with wrap around situation ---
         if (nLB == 12) then 
            if (tLB <= t ) then
               tUB = tUB + daysPerYear ! put UB in year 2
            else if (t < tUB ) then
               tLB = tLB - daysPerYear ! put LB in year 1
            else 
               call endrun(subName//"ERROR: in case A")
            end if
         end if
         !--- prevent infinite search ---
         n = n + 1
         if (n > 12) then
             write(iulog,F01) "ERROR: date,tLB,t,tUB = ",kyr,kmo,kda,tLB,t,tUB
             call endrun(subName//"ERROR: loop over first year, fail to find LB,UB")
         end if
      end do
   else if ( case == "C" ) then 
      !--- CASE C: loop over last year of data -----------------------
      if ( firstCallC ) then
         if ( masterproc ) write(iulog,F00) "CASE C: loop over last year of data"
         if ( masterproc ) write(iulog,F01) "CASE C: model date, t, time(nt) = ",kyr,kmo,kda,t,time(nt)
         nLB = nt-12 ; tLB = -2.0
         nUB = nt-11 ; tUB = -1.0 ! forces search for new LB,UB
         firstCallC = .false.
      end if
      t = mod(t,daysPerYear) + daysPerYear ! set year to 1
      n = 0
      readNewData = .false.
      do while (t < tLB  .or.  tUB < t)
         readNewData = .true.
         !--- move tUB,tLB forward in time ---
         nLB = nLB + 1 ; if (nLB > nt) nLB = nt - 11
         nUB = nLB + 1 ; if (nUB > nt) nUB = nt - 11
         tLB = mod(time(nLB),daysPerYear) + daysPerYear ! set year to 1
         tUB = mod(time(nUB),daysPerYear) + daysPerYear 
         !--- deal with wrap around situation ---
         if (nLB == nt) then 
            if (tLB <= t ) then
               tUB = tUB + daysPerYear ! put UB in year 2
               else if (t < tUB ) then
               tLB = tLB - daysPerYear ! put LB in year 1
            else 
               call endrun(subName//"ERROR: in case A")
            end if
         end if
         !--- prevent infinite search ---
         n = n + 1
         if (n > 12) then
             if ( masterproc ) write(iulog,F01) "ERROR: date,tLB,t,tUB = ",kyr,kmo,kda,tLB,t,tUB
             call endrun(subName//"ERROR: loop over first year, fail to find LB,UB")
         end if
      end do
   else
      !--- CASE B: interpolate within time series --------------------
      if ( firstCallB ) then
         if ( masterproc ) write(iulog,F01) "CASE B: interpolate within time series"
         if ( masterproc ) write(iulog,F01) "CASE B: date, time(1), model t, time(nt) = ",kyr,kmo,kda,time(1),t,time(nt)
         nLB = 0 ; tLB = -2.0
         nUB = 1 ; tUB = -1.0 ! forces search for new LB,UB
         firstCallB = .false.
      end if
      readNewData = .false.
      do while (tUB < t) 
         readNewData = .true.
         nLB = nLB + 1
         nUB = nLB + 1
         tLB = time(nLB)
         tUB = time(nUB)
         if (nUB > nt) call endrun(subName//"ERROR: nt < nUB")
      end do
   end if

   if (readNewData) then
      if ( masterproc ) write(iulog,F02) "read new data: model date, t = ", kyr,kmo,kda,ksec,t
      call readMonthlyAerdep (faerdep, nLB, nUB) ! input the new LB,UB data
   end if

   !----------------------------------------------------------------------------
   ! interpolate aerosol deposition data into 'forcing' array:
   !----------------------------------------------------------------------------
   fLB = (tUB - t)/(tUB - tLB)
   fUB = 1.0_r8 - fLB
   if (debug>2 .or. (debug==1 .and. (kda==1 .or. kda==15) .and. ksec==0)) then
      if ( masterproc ) write(iulog,F01) "date,tLB,t,tUB,nLB,nUB,fLB,fUB = ", kyr,kmo,kda,tLB,t,tUB,nLB,nUB,fLB,fUB
   endif

    do g = begg, endg
       if ( set_caerdep_from_file )then
          forc_aer(g, 1) = fLB*bcphidry2t(g,1)  + fUB*bcphidry2t(g,2)
          forc_aer(g, 2) = fLB*bcphodry2t(g,1)  + fUB*bcphodry2t(g,2)
          forc_aer(g, 3) = fLB*bcphiwet2t(g,1)  + fUB*bcphiwet2t(g,2)
          forc_aer(g, 4) = fLB*ocphidry2t(g,1)  + fUB*ocphidry2t(g,2)
          forc_aer(g, 5) = fLB*ocphodry2t(g,1)  + fUB*ocphodry2t(g,2)
          forc_aer(g, 6) = fLB*ocphiwet2t(g,1)  + fUB*ocphiwet2t(g,2)
       end if
       if ( set_dustdep_from_file )then
          forc_aer(g, 7) = fLB*dstx01wd2t(g,1)  + fUB*dstx01wd2t(g,2)
          forc_aer(g, 8) = fLB*dstx01dd2t(g,1)  + fUB*dstx01dd2t(g,2)
          forc_aer(g, 9) = fLB*dstx02wd2t(g,1)  + fUB*dstx02wd2t(g,2)
          forc_aer(g,10) = fLB*dstx02dd2t(g,1)  + fUB*dstx02dd2t(g,2)
          forc_aer(g,11) = fLB*dstx03wd2t(g,1)  + fUB*dstx03wd2t(g,2)
          forc_aer(g,12) = fLB*dstx03dd2t(g,1)  + fUB*dstx03dd2t(g,2)
          forc_aer(g,13) = fLB*dstx04wd2t(g,1)  + fUB*dstx04wd2t(g,2)
          forc_aer(g,14) = fLB*dstx04dd2t(g,1)  + fUB*dstx04dd2t(g,2)
       end if
    enddo


  end subroutine interpMonthlyAerdep

!================================================================================
!BOP
!
! !IROUTINE: readMonthlyAerdep
!
! !INTERFACE:
  subroutine readMonthlyAerdep (faer, n1, n2 )
!
! !DESCRIPTION:
! Read monthly aerosol deposition data for two consec. months.
!
! !USES:
    use clmtype
    use decompMod   , only : get_proc_bounds
    use clm_varpar  , only : lsmlon, lsmlat
    use clm_time_manager, only : get_nstep
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    character(len=*), intent(in) :: faer       ! file with monthly aerosol deposition
    integer, intent(in)          :: n1,n2      ! time index for samples 1 & 2
!
! !REVISION HISTORY:
! Adapted by Mark Flanner
!
!
!
! local pointers to implicit out arguments
!


! !LOCAL VARIABLES:
!EOP
    character(len=256) :: locfn           ! local file name
    integer :: k,n                        ! indices
    integer :: ncid,dimid,varid           ! input netCDF id's
    integer :: beg3d(3),len3d(3)          ! netCDF variable edges
    integer :: ntim                       ! number of input data time samples
    integer :: nlon_i                     ! number of input data longitudes
    integer :: nlat_i                     ! number of input data latitudes
    integer :: begg,endg                  ! beg and end local g index
    integer :: ier,ret                    ! error code
    integer :: closelatidx,closelonidx
    real(r8):: closelat,closelon

    real(r8), pointer :: arrayl(:)  ! temp local array

    character(*),parameter :: subname = '(readMonthlyAerdep) '
    
!--------------------------------------------------------------------------------

    ! Determine necessary indices

    call get_proc_bounds(begg=begg,endg=endg)

    ! ----------------------------------------------------------------------
    ! Open monthly aerosol file
    ! Read data and store in clmtype
    ! ----------------------------------------------------------------------

    if ( masterproc ) write(iulog,*) subName,'read samples ',n1,n2,' from ',trim(faer)

    do k=1,2   !loop over months and read aerosol data

       if (k==1) n = n1
       if (k==2) n = n2
       if (masterproc) then

          call getfil(faer, locfn, 0)
          call check_ret(nf_open(locfn, 0, ncid), subname)

          call check_ret(nf_inq_dimid (ncid, 'lon', dimid), subname)
          call check_ret(nf_inq_dimlen(ncid, dimid, nlon_i), subname)
          if (.not.single_column) then
             if (nlon_i /= lsmlon) then
                if ( masterproc ) write(iulog,*)subname,' parameter lsmlon= ',lsmlon,'does not equal input nlon_i= ',nlon_i
                call endrun( subname//' ERROR:: lsmlon does NOT equal input nlon' )
             end if
          endif
          call check_ret(nf_inq_dimid(ncid, 'lat', dimid), subname)
          call check_ret(nf_inq_dimlen(ncid, dimid, nlat_i), subname)

          if (.not.single_column ) then
             if (nlat_i /= lsmlat) then
                if ( masterproc ) write(iulog,*)subname,' parameter lsmlat= ',lsmlat,'does not equal input nlat_i= ',nlat_i
                call endrun( subname//' ERROR:: lsmlat does NOT equal input nlat' )
             end if
          endif

          call check_ret(nf_inq_dimid(ncid, 'time', dimid), subname)
          call check_ret(nf_inq_dimlen(ncid, dimid, ntim), subname)

       else

          nlon_i = lsmlon
          nlat_i = lsmlon
 
       endif   ! masterproc

       call mpi_bcast (nlon_i, 1, MPI_INTEGER, 0, mpicom, ier)
       call mpi_bcast (nlat_i, 1, MPI_INTEGER, 0, mpicom, ier)
       call mpi_bcast (ntim  , 1, MPI_INTEGER, 0, mpicom, ier)

       if (single_column) then
          call scam_setlatlonidx(ncid,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
       endif

       allocate(arrayl(begg:endg),stat=ier)
       if (ier /= 0) call endrun( subname//' ERROR:: allocation array error' )

       if (single_column) then
          beg3d(1) = closelonidx; len3d(1) = 1
          beg3d(2) = closelatidx; len3d(2) = 1
       else
          beg3d(1) = 1         ; len3d(1) = nlon_i
          beg3d(2) = 1         ; len3d(2) = nlat_i
       end if

       beg3d(3) = n ; len3d(3) = 1

       
       call ncd_iolocal(ncid,'BCDEPWET','read',arrayl,grlnd,beg3d,len3d,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: MONTHLY_BCPHIWET NOT on faer file' )
       bcphiwet2t(begg:endg,k) = arrayl(begg:endg)

       call ncd_iolocal(ncid,'BCPHIDRY','read',arrayl,grlnd,beg3d,len3d,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: MONTHLY_BCPHIDRY NOT on faer file' )
       bcphidry2t(begg:endg,k) = arrayl(begg:endg)

       call ncd_iolocal(ncid,'BCPHODRY','read',arrayl,grlnd,beg3d,len3d,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: MONTHLY_BCPHODRY NOT on faer file' )
       bcphodry2t(begg:endg,k) = arrayl(begg:endg)

       call ncd_iolocal(ncid,'OCDEPWET','read',arrayl,grlnd,beg3d,len3d,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: MONTHLY_OCPHIWET NOT on faer file' )
       ocphiwet2t(begg:endg,k) = arrayl(begg:endg)

       call ncd_iolocal(ncid,'OCPHIDRY','read',arrayl,grlnd,beg3d,len3d,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: MONTHLY_OCPHIDRY NOT on faer file' )
       ocphidry2t(begg:endg,k) = arrayl(begg:endg)

       call ncd_iolocal(ncid,'OCPHODRY','read',arrayl,grlnd,beg3d,len3d,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: MONTHLY_OCPHODRY NOT on faer file' )
       ocphodry2t(begg:endg,k) = arrayl(begg:endg)

       call ncd_iolocal(ncid,'DSTX01WD','read',arrayl,grlnd,beg3d,len3d,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: MONTHLY_DSTX01WD NOT on faer file' )
       dstx01wd2t(begg:endg,k) = arrayl(begg:endg)

       call ncd_iolocal(ncid,'DSTX01DD','read',arrayl,grlnd,beg3d,len3d,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: MONTHLY_DSTX01DD NOT on faer file' )
       dstx01dd2t(begg:endg,k) = arrayl(begg:endg)

       call ncd_iolocal(ncid,'DSTX02WD','read',arrayl,grlnd,beg3d,len3d,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: MONTHLY_DSTX02WD NOT on faer file' )
       dstx02wd2t(begg:endg,k) = arrayl(begg:endg)

       call ncd_iolocal(ncid,'DSTX02DD','read',arrayl,grlnd,beg3d,len3d,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: MONTHLY_DSTX02DD NOT on faer file' )
       dstx02dd2t(begg:endg,k) = arrayl(begg:endg)
       
       call ncd_iolocal(ncid,'DSTX03WD','read',arrayl,grlnd,beg3d,len3d,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: MONTHLY_DSTX03WD NOT on faer file' )
       dstx03wd2t(begg:endg,k) = arrayl(begg:endg)

       call ncd_iolocal(ncid,'DSTX03DD','read',arrayl,grlnd,beg3d,len3d,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: MONTHLY_DSTX03DD NOT on faer file' )
       dstx03dd2t(begg:endg,k) = arrayl(begg:endg)

       call ncd_iolocal(ncid,'DSTX04WD','read',arrayl,grlnd,beg3d,len3d,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: MONTHLY_DSTX04WD NOT on faer file' )
       dstx04wd2t(begg:endg,k) = arrayl(begg:endg)

       call ncd_iolocal(ncid,'DSTX04DD','read',arrayl,grlnd,beg3d,len3d,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: MONTHLY_DSTX04DD NOT on faer file' )
       dstx04dd2t(begg:endg,k) = arrayl(begg:endg)

       deallocate(arrayl)

       if (masterproc) then
          call check_ret(nf_close(ncid), subname)
          if (debug>1) then
             write(iulog,*) subName,'Successfully read aerosol data for index n = ',n
             write(iulog,*)
             call shr_sys_flush(iulog)
          end if
       end if

    end do   ! end of loop over months

  end subroutine readMonthlyAerdep

end module aerdepMod
