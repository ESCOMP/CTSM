#include <misc.h>
#include <preproc.h>
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
  use perf_mod,        only : t_startf, t_stopf
  use clm_varctl,      only : set_caerdep_from_file, set_dustdep_from_file
!
! !PUBLIC TYPES:
  implicit none
  save
 
  private
  
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: interpMonthlyAerdep   ! interpolate monthly aerosol deposition data
  public :: aerdepini             ! aerosol deposition initialization

!
! !REVISION HISTORY:
! Created by Mark Flanner, 
!   based on vegetation interpolation schemes in STATICEcosystemDynMod
!
!EOP
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: readMonthlyAerdep       ! read monthly aerosol deposition data for two months

!
! PRIVATE TYPES:
  integer , private :: InterpMonths1_aer        ! saved month index (for aerosol deposition)

  real(r8), private, allocatable :: bcphiwet2t(:,:)
  real(r8), private, allocatable :: bcphidry2t(:,:)
  real(r8), private, allocatable :: bcphodry2t(:,:)
  real(r8), private, allocatable :: ocphiwet2t(:,:)
  real(r8), private, allocatable :: ocphidry2t(:,:)
  real(r8), private, allocatable :: ocphodry2t(:,:)
  real(r8), private, allocatable :: dstx01wd2t(:,:)
  real(r8), private, allocatable :: dstx01dd2t(:,:)
  real(r8), private, allocatable :: dstx02wd2t(:,:)
  real(r8), private, allocatable :: dstx02dd2t(:,:)
  real(r8), private, allocatable :: dstx03wd2t(:,:)
  real(r8), private, allocatable :: dstx03dd2t(:,:)
  real(r8), private, allocatable :: dstx04wd2t(:,:)
  real(r8), private, allocatable :: dstx04dd2t(:,:)

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
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
!
! !ARGUMENTS:
    implicit none

!
! !REVISION HISTORY:
!
!EOP
!
! LOCAL VARIABLES:
    integer :: ier         ! error code
    integer :: begg,endg   ! local beg and end p index
!-----------------------------------------------------------------------

    InterpMonths1_aer = -999  ! saved month index

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
          
       if (ier /= 0) then
          write(iulog,*) 'aerdepini allocation error'
          call endrun
       end if

        bcphiwet2t(begg:endg,1:2) = nan
        bcphidry2t(begg:endg,1:2) = nan
        bcphodry2t(begg:endg,1:2) = nan
        ocphiwet2t(begg:endg,1:2) = nan
        ocphidry2t(begg:endg,1:2) = nan
        ocphodry2t(begg:endg,1:2) = nan
    end if

    if ( set_dustdep_from_file )then

       allocate(dstx01wd2t(begg:endg,2))
       allocate(dstx01dd2t(begg:endg,2))
       allocate(dstx02wd2t(begg:endg,2))
       allocate(dstx02dd2t(begg:endg,2))
       allocate(dstx03wd2t(begg:endg,2))
       allocate(dstx03dd2t(begg:endg,2))
       allocate(dstx04wd2t(begg:endg,2))
       allocate(dstx04dd2t(begg:endg,2))

       if (ier /= 0) then
          write(iulog,*) 'aerdepini allocation error'
          call endrun
       end if

       dstx01wd2t(begg:endg,1:2) = nan
       dstx01dd2t(begg:endg,1:2) = nan
       dstx02wd2t(begg:endg,1:2) = nan
       dstx02dd2t(begg:endg,1:2) = nan
       dstx03wd2t(begg:endg,1:2) = nan
       dstx03dd2t(begg:endg,1:2) = nan
       dstx04wd2t(begg:endg,1:2) = nan
       dstx04dd2t(begg:endg,1:2) = nan

    end if

  end subroutine aerdepini


!-----------------------------------------------------------------------
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
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!  Adapted by Mark Flanner
!
!EOP
!

!
! local pointers to implicit out arguments
!
    real(r8), pointer :: forc_aer(:,:)   ! aerosol deposition rate (kg/m2/s)


! LOCAL VARIABLES:
    real(r8):: timwt_aer(2)  ! time weights for month 1 and month 2 (aerosol deposition)
    integer :: kyr         ! year (0, ...) for nstep+1
    integer :: kmo         ! month (1, ..., 12)
    integer :: kda         ! day of month (1, ..., 31)
    integer :: ksec        ! seconds into current date for nstep+1
    real(r8):: dtime       ! land model time step (sec)
    real(r8):: t           ! a fraction: kda/ndaypm
    integer :: it(2)       ! month 1 and month 2 (step 1)
    integer :: months(2)   ! months to be interpolated (1 to 12)
    integer :: g
    integer :: begg,endg                  ! beg and end local g index
    
    integer, dimension(12) :: ndaypm= &
         (/31,28,31,30,31,30,31,31,30,31,30,31/) !days per month
!-----------------------------------------------------------------------


    ! Determine necessary indices
    call get_proc_bounds(begg=begg,endg=endg)
    

    ! Assign local pointers to derived subtypes components (gridcell level)
    forc_aer     => clm_a2l%forc_aer
    
   
    ! get timestep
    dtime = get_step_size()

    if ( is_perpetual() ) then
       call get_perp_date(kyr, kmo, kda, ksec, offset=int(dtime))
    else
       call get_curr_date(kyr, kmo, kda, ksec, offset=int(dtime))
    end if

    t = (kda-0.5_r8) / ndaypm(kmo)
    it(1) = t + 0.5_r8
    it(2) = it(1) + 1
    months(1) = kmo + it(1) - 1
    months(2) = kmo + it(2) - 1
    if (months(1) <  1) months(1) = 12
    if (months(2) > 12) months(2) = 1
    timwt_aer(1) = (it(1)+0.5_r8) - t
    timwt_aer(2) = 1._r8-timwt_aer(1)

    if (InterpMonths1_aer /= months(1)) then
       call t_startf('readMonthlyAerdep')
       call readMonthlyAerdep (faerdep, kmo, kda, months)
       InterpMonths1_aer = months(1)
       call t_stopf('readMonthlyAerdep')
    end if


    ! interpolate aerosol deposition data into 'forcing' array:
    do g = begg, endg
       if ( set_caerdep_from_file )then
          forc_aer(g,1) = timwt_aer(1)*bcphidry2t(g,1)  + timwt_aer(2)*bcphidry2t(g,2)
          forc_aer(g,2) = timwt_aer(1)*bcphodry2t(g,1)  + timwt_aer(2)*bcphodry2t(g,2)
          forc_aer(g,3) = timwt_aer(1)*bcphiwet2t(g,1)  + timwt_aer(2)*bcphiwet2t(g,2)
          forc_aer(g,4) = timwt_aer(1)*ocphidry2t(g,1)  + timwt_aer(2)*ocphidry2t(g,2)
          forc_aer(g,5) = timwt_aer(1)*ocphodry2t(g,1)  + timwt_aer(2)*ocphodry2t(g,2)
          forc_aer(g,6) = timwt_aer(1)*ocphiwet2t(g,1)  + timwt_aer(2)*ocphiwet2t(g,2)
       end if
       if ( set_dustdep_from_file )then
          forc_aer(g,7) = timwt_aer(1)*dstx01wd2t(g,1)  + timwt_aer(2)*dstx01wd2t(g,2)
          forc_aer(g,8) = timwt_aer(1)*dstx01dd2t(g,1)  + timwt_aer(2)*dstx01dd2t(g,2)
          forc_aer(g,9) = timwt_aer(1)*dstx02wd2t(g,1)  + timwt_aer(2)*dstx02wd2t(g,2)
          forc_aer(g,10) = timwt_aer(1)*dstx02dd2t(g,1)  + timwt_aer(2)*dstx02dd2t(g,2)
          forc_aer(g,11) = timwt_aer(1)*dstx03wd2t(g,1)  + timwt_aer(2)*dstx03wd2t(g,2)
          forc_aer(g,12) = timwt_aer(1)*dstx03dd2t(g,1)  + timwt_aer(2)*dstx03dd2t(g,2)
          forc_aer(g,13) = timwt_aer(1)*dstx04wd2t(g,1)  + timwt_aer(2)*dstx04wd2t(g,2)
          forc_aer(g,14) = timwt_aer(1)*dstx04dd2t(g,1)  + timwt_aer(2)*dstx04dd2t(g,2)
       end if
    enddo

  end subroutine interpMonthlyAerdep


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readMonthlyAerdep
!
! !INTERFACE:
  subroutine readMonthlyAerdep (faer, kmo, kda, months)
!
! !DESCRIPTION:
! Read monthly aerosol deposition data for two consec. months.
!
! !USES:
    use clmtype
    use decompMod   , only : get_proc_bounds
    use clm_varpar  , only : lsmlon, lsmlat
    use fileutils   , only : getfil
    use spmdMod     , only : masterproc, mpicom, MPI_REAL8, MPI_INTEGER
    use clm_time_manager, only : get_nstep
    use ncdio       , only : check_ret,ncd_iolocal
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    character(len=*), intent(in) :: faer       ! file with monthly aerosol deposition
    integer, intent(in)          :: kmo        ! month (1, ..., 12)
    integer, intent(in)          :: kda        ! day of month (1, ..., 31)
    integer, intent(in)          :: months(2)  ! months to be interpolated (1 to 12)
!
! !REVISION HISTORY:
! Adapted by Mark Flanner
!
!EOP
!
!
! local pointers to implicit out arguments
!


! LOCAL VARIABLES:
    character(len=256) :: locfn           ! local file name
    integer :: k                          ! indices
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

    character(len=32) :: subname = 'readMonthlyAerdep'

    
!-----------------------------------------------------------------------

    ! Determine necessary indices

    call get_proc_bounds(begg=begg,endg=endg)

    ! ----------------------------------------------------------------------
    ! Open monthly aerosol file
    ! Read data and store in clmtype
    ! ----------------------------------------------------------------------

    do k=1,2   !loop over months and read aerosol data

       if (masterproc) then

          write(iulog,*) 'Attempting to read monthly aerosol deposition data .....'
          write(iulog,*) 'nstep = ',get_nstep(),' month = ',kmo,' day = ',kda

          call getfil(faer, locfn, 0)
          call check_ret(nf_open(locfn, 0, ncid), subname)

          call check_ret(nf_inq_dimid (ncid, 'lon', dimid), subname)
          call check_ret(nf_inq_dimlen(ncid, dimid, nlon_i), subname)
          if (.not.single_column) then
             if (nlon_i /= lsmlon) then
                write(iulog,*)subname,' parameter lsmlon= ',lsmlon,'does not equal input nlon_i= ',nlon_i
                call endrun()
             end if
          endif
          call check_ret(nf_inq_dimid(ncid, 'lat', dimid), subname)
          call check_ret(nf_inq_dimlen(ncid, dimid, nlat_i), subname)

          if (.not.single_column ) then
             if (nlat_i /= lsmlat) then
                write(iulog,*)subname,' parameter lsmlat= ',lsmlat,'does not equal input nlat_i= ',nlat_i
                call endrun()
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
       if (ier /= 0) then
          write(iulog,*)subname, 'allocation array error '; call endrun()
       end if

       if (single_column) then
          beg3d(1) = closelonidx; len3d(1) = 1
          beg3d(2) = closelatidx; len3d(2) = 1
       else
          beg3d(1) = 1         ; len3d(1) = nlon_i
          beg3d(2) = 1         ; len3d(2) = nlat_i
       end if

       beg3d(3) = months(k) ; len3d(3) = 1

       
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
          write(iulog,*) 'Successfully read aerosol deposition data for'
          write(iulog,*) 'month ', months(k)
          write(iulog,*)
       end if

    end do   ! end of loop over months

  end subroutine readMonthlyAerdep

end module aerdepMod
